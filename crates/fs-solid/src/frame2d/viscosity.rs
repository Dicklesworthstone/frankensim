//! Original material section viscosity on the native frame/Newmark trajectory.
//!
//! Integrate axial strain rate and cubic-Hermite curvature rate in the same
//! reference member. The full assembled C retains joint and directional cross
//! terms; heterogeneous members need not share one proportional damping ratio.
//! Material work is exported as dissipation, not silently deposited into heat.

use super::{Frame2d, Geometry, SectionModel, SolidError};
use fs_material::section::{
    ConcentricElasticSection, ConcentricKelvinVoigtSection, axial::ConcentricAxialViscosity,
};

/// Explicit axial AND bending viscosities from one original material section.
/// Each original layer/source and its frequency applicability remain attached.
#[derive(Clone, Debug)]
pub struct FrameViscousSection {
    axial: ConcentricAxialViscosity,
    bending: ConcentricKelvinVoigtSection,
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
    maximum_omega_rad_s: f64,
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
    /// axial, symmetric bending, antisymmetric bending. Explicit zero factors
    /// are retained, and `C = sum(damping_n_s_m * weights * weights^T)`.
    pub fn ports(&self) -> &[FrameViscousPort] {
        &self.ports
    }
    /// Original unconstrained full-coordinate indices; no condensation is used.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Original member index and both original layered source bindings.
    pub fn sections(&self) -> &[(usize, FrameViscousSection)] {
        &self.sections
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
        if self.maximum_omega_rad_s >= nyquist {
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
            band,
        })
    }

    /// Same geometry, elastic stiffness and inertia used by both loss owners.
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
    /// Intersection of both original coefficient bands [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.band
    }
}

impl Frame2d {
    /// Add an elastic member whose SAME sources also determine its physical
    /// strain-rate damping. Static equilibrium and mass are unchanged by eta;
    /// every later native Newmark call automatically includes these losses.
    /// The elastic coefficients retain `add_material_element`'s unrestricted
    /// frequency support, so static and ordinary undamped modal queries remain
    /// valid. Both explicit viscosity laws retain their separate finite bands.
    ///
    /// This is frozen-temperature, small-strain linear geometry, with no
    /// Maxwell history, creep, thermal feedback or radiation model. All frame
    /// members must remain linear elastic during a material-damped history.
    /// Free coordinates must have positive reference inertia; choose consistent
    /// member mass when free rotations need inertia. Static massless condensation
    /// cannot represent a viscous first-order degree of freedom.
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
        let element = self.add_material_element(i, j, np, section.elastic().clone())?;
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
        let mut maximum_omega_rad_s = 0.0_f64;
        for mode in self.modal_analysis()?.modes {
            let omega = mode.angular_frequency().ok_or_else(|| {
                invalid("complete material-damped frame spectrum must be positive and resolved")
            })?;
            maximum_omega_rad_s = maximum_omega_rad_s.max(omega);
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
            let (l, c, s) = (element.l0, element.cos0, element.sin0);
            let ca = section.axial.viscous_rigidity_n_s() / l;
            let cb = section.bending.viscous_stiffness_n_m2_s() / l / l / l;
            if (section.axial.viscous_rigidity_n_s() > 0.0 && ca == 0.0)
                || (section.bending.viscous_stiffness_n_m2_s() > 0.0 && cb == 0.0)
                || !ca.is_finite()
                || !(3.0 * cb).is_finite()
            {
                return Err(invalid(
                    "material section viscosity is unrepresentable on the member geometry",
                ));
            }
            let dofs = [
                3 * element.i,
                3 * element.i + 1,
                3 * element.i + 2,
                3 * element.j,
                3 * element.j + 1,
                3 * element.j + 2,
            ];
            for (damping_n_s_m, displacement_weights) in [
                (ca, [-c, -s, 0.0, c, s, 0.0]),
                (3.0 * cb, [-2.0 * s, 2.0 * c, l, 2.0 * s, -2.0 * c, l]),
                (cb, [0.0, 0.0, l, 0.0, 0.0, -l]),
            ] {
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
            maximum_omega_rad_s,
        }))
    }
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
