//! Sparse signed spectra about the SAME loaded finite-strain sheet equilibrium.
//!
//! This adapter keeps the native plane-stress/exact-volume or mean-dilatation
//! tangent and original consistent reference mass sparse. The existing fs-modal
//! owner supplies shared sparse factors, inertia counts and complete-in-window
//! residual certificates; no dense n-by-n matrix or 128-free-DOF cap is added.
//!
//! Windows use signed squared-frequency units [s^-2] and the native (low, high]
//! convention. Negative modes remain negative. Positive selected windows do
//! not establish stability outside their coverage or a dynamics truncation.
//! Singular endpoints refuse without being nudged. The native residual/|lambda|
//! gate can refuse nearly neutral/rigid modes; the complete bounded dense route
//! remains available for small pencils. These remain in-plane modes with
//! quasistatic transverse relaxation, not bending or wrinkling modes.
//!
//! The explicit mass budget bounds original inertia assembly/retention. The
//! existing Newton/constitutive owners retain equilibrium and tangent work;
//! WindowPlan and SliceOptions separately bound the free-pencil spectral
//! payload, sparse factors, iterations and workers. They do not bound total
//! process memory or certify a particular problem size/performance.

use super::{
    PlaneReferenceDensity, PlaneVibrationError, PlaneVibrationFormulation, PlaneVibrationState,
    checkpoint, inertia_error, invalid,
};
use crate::hyper2d::{HyperProblem, NewtonReport};
use fs_exec::Cx;
use fs_modal::{MultiSliceReport, SliceOptions, WindowPlan, slice_windows};
use fs_sparse::Csr;

/// Original loaded equilibrium, reference inertia and every mode in the
/// requested signed windows. No claims are made for omitted eigenvalues.
pub struct PlaneVibrationWindows<'p, 'm> {
    state: PlaneVibrationState<'p, 'm>,
    stiffness_n_m: Csr,
    mass_kg: Csr,
    plan: WindowPlan,
    options: SliceOptions,
    spectrum: MultiSliceReport,
}

// The parent prepares canonical native rows and an increasing free map.
// Restrict their sparse rows directly, without COO sorting,
// dense probing, mass lumping or a separately reconstructed prestress matrix.
fn restrict(
    matrix: &Csr,
    free: &[usize],
    reduced: &[usize],
    scale: f64,
    cx: &Cx<'_>,
) -> Result<Csr, crate::SolidError> {
    let n = free.len();
    let mut rows = Vec::with_capacity(n + 1);
    let mut columns = Vec::new();
    let mut values = Vec::new();
    let mut largest = 0.0_f64;
    rows.push(0);
    for &row in free {
        checkpoint(cx)?;
        let (native_columns, native_values) = matrix.row(row);
        for (&column, &value) in native_columns.iter().zip(native_values) {
            let column = reduced[column];
            if column == usize::MAX {
                continue;
            }
            let scaled = value * scale;
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                return Err(invalid("unrepresentable physical plane vibration pencil"));
            }
            columns.push(column);
            values.push(scaled);
            largest = largest.max(scaled.abs());
        }
        rows.push(values.len());
    }
    let result =
        Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, || checkpoint(cx))?
            .ok_or_else(|| invalid("noncanonical native plane vibration restriction"))?;
    let tolerance = (64.0 * f64::EPSILON * n as f64) * largest;
    for row in 0..n {
        checkpoint(cx)?;
        let (columns, values) = result.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            if (value - result.get(column, row)).abs() > tolerance {
                return Err(invalid(
                    "plane vibration pencil is not symmetric to roundoff",
                ));
            }
        }
    }
    Ok(result)
}

impl<'m> HyperProblem<'m> {
    /// Solve the original loaded equilibrium and every requested sparse window.
    ///
    /// Original density [kg/m³] and thickness [m] retain their source frame and
    /// mass inventory. The selected native owner supplies equilibrium and its
    /// complete condensed Hessian, including initial stress, exactly once.
    /// The mass-assembly byte limit is independent of the spectral budgets.
    ///
    /// Optional geometric-ordering coordinates in the plan use ascending free
    /// node-major xy order: repeat [X,Y,0] for each free component at that node.
    /// Native material, reference, support, mass, symmetry, spectral budget,
    /// unresolved-window and cancellation refusals return no partial result.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_vibration_windows<'p>(
        &'p self,
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        maximum_mass_assembly_bytes: usize,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibrationWindows<'p, 'm>, PlaneVibrationError> {
        self.solve_vibration_windows_impl(
            formulation,
            reference_thickness_m,
            density,
            None,
            maximum_mass_assembly_bytes,
            plan,
            options,
            cx,
        )
    }

    /// Same sparse physical solve with a Newton starting displacement.
    /// The original absolute support targets and loads are re-equilibrated;
    /// the starting field never bypasses equilibrium or supplies a foreign K.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_vibration_windows_from<'p>(
        &'p self,
        initial_displacement: &[[f64; 2]],
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        maximum_mass_assembly_bytes: usize,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibrationWindows<'p, 'm>, PlaneVibrationError> {
        self.solve_vibration_windows_impl(
            formulation,
            reference_thickness_m,
            density,
            Some(initial_displacement),
            maximum_mass_assembly_bytes,
            plan,
            options,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_vibration_windows_impl<'p>(
        &'p self,
        formulation: PlaneVibrationFormulation,
        reference_thickness_m: f64,
        density: PlaneReferenceDensity<'p>,
        initial: Option<&[[f64; 2]]>,
        maximum_mass_assembly_bytes: usize,
        plan: &WindowPlan,
        options: &SliceOptions,
        cx: &Cx<'_>,
    ) -> Result<PlaneVibrationWindows<'p, 'm>, PlaneVibrationError> {
        let (state, tangent) = self.equilibrate_for_vibration(
            formulation,
            reference_thickness_m,
            density,
            initial,
            usize::MAX,
            maximum_mass_assembly_bytes,
            cx,
        )?;
        let mut reduced = vec![usize::MAX; 2 * self.mesh.nodes.len()];
        for (index, &dof) in state.free_dofs.iter().enumerate() {
            reduced[dof] = index;
        }
        let stiffness_n_m = restrict(
            &tangent,
            &state.free_dofs,
            &reduced,
            reference_thickness_m,
            cx,
        )?;
        let mass_kg = restrict(
            state.reference_mass.matrix(),
            &state.free_dofs,
            &reduced,
            1.0,
            cx,
        )?;
        // Release full tangent and restriction scratch before spectral admission.
        // The full reference mass remains owned for physical acceleration loads.
        drop(tangent);
        drop(reduced);
        let spectrum = slice_windows(&stiffness_n_m, &mass_kg, plan, options, cx)
            .map_err(PlaneVibrationError::Modal)?;
        let result = PlaneVibrationWindows {
            state,
            stiffness_n_m,
            mass_kg,
            plan: plan.clone(),
            options: *options,
            spectrum,
        };
        checkpoint(cx)?;
        Ok(result)
    }
}

impl<'p, 'm> PlaneVibrationWindows<'p, 'm> {
    /// Original geometry, material/reference owners, supports and dead loads.
    pub fn problem(&self) -> &HyperProblem<'m> {
        self.state.problem
    }
    pub const fn formulation(&self) -> PlaneVibrationFormulation {
        self.state.formulation
    }
    pub const fn reference_thickness_m(&self) -> f64 {
        self.state.reference_thickness_m
    }
    pub const fn reference_density(&self) -> PlaneReferenceDensity<'p> {
        self.state.density
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.state.reference_mass.reference_mass_kg()
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.state.displacement_m
    }
    pub fn equilibrium_report(&self) -> &NewtonReport {
        &self.state.equilibrium
    }
    /// Signed forces BY the supports ON the sheet in the original nodal frame.
    pub fn reactions_n(&self) -> &[(usize, f64)] {
        &self.state.reactions_n
    }
    /// Ascending reduced-to-full node-major xy DOFs; prescribed increments vanish.
    pub fn free_dofs(&self) -> &[usize] {
        &self.state.free_dofs
    }
    /// Original full condensed free-coordinate tangent, including initial stress.
    pub const fn stiffness_n_m(&self) -> &Csr {
        &self.stiffness_n_m
    }
    /// Original free-coordinate consistent inertia, without density averaging.
    pub const fn mass_kg(&self) -> &Csr {
        &self.mass_kg
    }
    pub const fn plan(&self) -> &WindowPlan {
        &self.plan
    }
    pub const fn options(&self) -> &SliceOptions {
        &self.options
    }
    /// Signed modes, endpoint inertia counts and complete-in-window certificates.
    pub const fn spectrum(&self) -> &MultiSliceReport {
        &self.spectrum
    }

    /// SAME full consistent mass action, including prescribed coordinates.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 2]],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, PlaneVibrationError> {
        self.state
            .reference_mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)
            .map_err(Into::into)
    }
}
