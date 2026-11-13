//! One physical mesh path; every sample uses the original dielectric FE solver.

use fs_em::electrostatic::{ElectrostaticControls, PlanarElectrostatics};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::field_storage::planar::PlanarDielectricLaw;
use fs_math::dd::Dd;
use fs_rep_mesh::{Metric2, TriComplex2};

use super::{Error, Result, check, finite, positive};

/// An affine nodal path `x(s) = x_reference + s * motion`, with s in metres.
///
/// The positive electrode translates rigidly in the unit `axis` direction;
/// the negative electrode is fixed. Interior mesh motion and any remaining
/// boundary motion are explicitly authored parts of the same physical path.
/// The caller must represent their mechanical work in this ONE coordinate;
/// this is not an automatic mesh-motion or collision algorithm.
///
/// Original isotropic dielectric states remain constant on their cells. The
/// constitutive assumption is strain-independent epsilon, with no density,
/// electrostriction, thermal evolution or dielectric loss contribution. Fixed
/// extrusion thickness omits end fields. Nonuniform planar fields, fringing
/// inside the supplied domain and interfaces are resolved by the original
/// Hodge/CG owner. Unlisted outer edges retain its insulating boundary.
pub struct PlanarCapacitorPath {
    reference: PlanarElectrostatics,
    motion: Vec<[f64; 2]>,
    terminals: [Vec<usize>; 2],
    limits_m: [f64; 2],
    axis: [f64; 2],
    assembly: WeightedAssemblyLimits,
    controls: ElectrostaticControls,
    reference_sample: Option<CapacitanceSample>,
}

/// Unit-voltage response from the actual field at one admitted geometry.
#[derive(Clone, Copy, Debug)]
pub struct CapacitanceSample {
    coordinate_m: f64,
    capacitance_f: f64,
    derivative_f_per_m: f64,
    charge_defect_per_v: f64,
    energy_defect_per_v2: f64,
}
impl CapacitanceSample {
    pub const fn coordinate_m(&self) -> f64 {
        self.coordinate_m
    }
    pub const fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
    /// Analytic shape derivative, evaluated on the same unit-voltage solution.
    pub const fn derivative_f_per_m(&self) -> f64 {
        self.derivative_f_per_m
    }
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

impl PlanarCapacitorPath {
    /// Bind one resolved field, original terminal vertex IDs and complete
    /// nodal motion. Open coordinate limits must straddle the reference s=0.
    /// Every cell's quadratic signed-area ratio is checked over the ENTIRE
    /// interval, including its interior minimum, before motion is admitted.
    /// This prevents local folding, not global self-intersection or contact.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: PlanarElectrostatics,
        terminals: [Vec<usize>; 2],
        motion: Vec<[f64; 2]>,
        limits_m: [f64; 2],
        assembly: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        let mesh = reference.mesh();
        if !matches!(mesh.metric(), Metric2::Planar { .. }) {
            return Err(Error::Invalid(
                "guided Cartesian capacitor motion requires a planar fixed-thickness field",
            ));
        }
        if motion.len() != mesh.vertices().len()
            || motion.iter().flatten().any(|x| !x.is_finite())
            || limits_m.iter().any(|x| !x.is_finite())
            || limits_m[0] >= 0.0
            || limits_m[1] <= 0.0
            || terminals.iter().any(Vec::is_empty)
            || terminals.iter().flatten().any(|&i| i >= motion.len())
            || mesh.vertices().len() > assembly.max_dofs
            || mesh.edges().len() > assembly.max_dofs
            || mesh.faces().len() > assembly.max_cells
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > assembly.max_triplets)
            || reference
                .permittivity_laws()
                .iter()
                .any(|law| !matches!(law, PlanarDielectricLaw::Isotropic))
        {
            return Err(Error::Invalid(
                "bounded complete isotropic mesh path, terminals and open limits required",
            ));
        }
        let axis = motion[terminals[0][0]];
        let length = finite(axis[0].hypot(axis[1]))?;
        if (length - 1.0).abs() > 16.0 * f64::EPSILON
            || terminals[0].iter().any(|&i| motion[i] != axis)
            || terminals[1].iter().any(|&i| motion[i] != [0.0; 2])
        {
            return Err(Error::Invalid(
                "one rigid unit-direction moving electrode and a fixed opposite electrode required",
            ));
        }
        for face in mesh.faces() {
            check(cancelled)?;
            let [a, b, c] = face.map(|i| mesh.vertices()[i as usize]);
            let [u, v, w] = face.map(|i| motion[i as usize]);
            let e1 = sub(b, a);
            let e2 = sub(c, a);
            let d1 = sub(v, u);
            let d2 = sub(w, u);
            let area = finite(cross(e1, e2))?;
            let linear = finite((cross(d1, e2) + cross(e1, d2)) / area)?;
            let quadratic = finite(cross(d1, d2) / area)?;
            let ratio = |s: f64| quadratic.mul_add(s, linear).mul_add(s, 1.0);
            let mut minimum = finite(ratio(limits_m[0]))?.min(finite(ratio(limits_m[1]))?);
            if quadratic > 0.0 {
                let vertex = -0.5 * (linear / quadratic);
                if vertex > limits_m[0] && vertex < limits_m[1] {
                    minimum = minimum.min(finite(ratio(vertex))?);
                }
            }
            if !positive(minimum) {
                return Err(Error::Invalid(
                    "electrostatic mesh path folds or loses area inside its coordinate interval",
                ));
            }
        }
        let mut path = Self {
            reference,
            motion,
            terminals,
            limits_m,
            axis,
            assembly,
            controls,
            reference_sample: None,
        };
        // The original owner admits topology, electrode coverage and gauges.
        path.reference_sample = Some(path.sample(0.0, cancelled)?);
        check(cancelled)?;
        Ok(path)
    }

    pub const fn reference(&self) -> &PlanarElectrostatics {
        &self.reference
    }
    pub fn nodal_motion_per_m(&self) -> &[[f64; 2]] {
        &self.motion
    }
    pub fn terminals(&self) -> &[Vec<usize>; 2] {
        &self.terminals
    }
    pub const fn coordinate_limits_m(&self) -> [f64; 2] {
        self.limits_m
    }
    pub const fn translation_axis(&self) -> [f64; 2] {
        self.axis
    }

    fn admit_coordinate(&self, coordinate: f64) -> Result<()> {
        if !coordinate.is_finite()
            || coordinate <= self.limits_m[0]
            || coordinate >= self.limits_m[1]
        {
            Err(Error::Invalid(
                "electrostatic coordinate is outside the admitted open mesh path",
            ))
        } else {
            Ok(())
        }
    }

    /// Rebuild the original source-bound operator at a physical coordinate.
    /// Vertex keys, lineage, connectivity, cells, laws and thickness survive.
    /// No interpolation of a capacitance table or fitted force law is used.
    pub fn problem_at(
        &self,
        coordinate: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarElectrostatics> {
        check(cancelled)?;
        self.admit_coordinate(coordinate)?;
        let old = self.reference.mesh();
        let mut positions = Vec::with_capacity(old.vertices().len());
        for (point, motion) in old.vertices().iter().zip(&self.motion) {
            check(cancelled)?;
            let mut next = [0.0; 2];
            for axis in 0..2 {
                next[axis] = finite(coordinate.mul_add(motion[axis], point[axis]))?;
                if coordinate != 0.0 && motion[axis] != 0.0 && next[axis] == point[axis] {
                    return Err(Error::Invalid(
                        "physical mesh motion disappeared in coordinate rounding",
                    ));
                }
            }
            positions.push(next);
        }
        let mesh = TriComplex2::from_triangles(
            old.lineage(),
            positions,
            old.vertex_keys().to_vec(),
            old.faces().to_vec(),
            old.metric(),
        )
        .map_err(|_| Error::Invalid("moved electrostatic mesh is not admissible"))?;
        Ok(PlanarElectrostatics::new_with_laws(
            mesh,
            self.reference.materials().to_vec(),
            self.reference.permittivity_laws().to_vec(),
            self.assembly,
            cancelled,
        )?)
    }

    pub fn sample(
        &self,
        coordinate: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CapacitanceSample> {
        check(cancelled)?;
        self.admit_coordinate(coordinate)?;
        if coordinate == 0.0 {
            match self.reference_sample {
                Some(sample) => Ok(sample),
                None => self.sample_problem(&self.reference, coordinate, cancelled),
            }
        } else {
            self.sample_problem(
                &self.problem_at(coordinate, cancelled)?,
                coordinate,
                cancelled,
            )
        }
    }

    fn sample_problem(
        &self,
        problem: &PlanarElectrostatics,
        coordinate: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CapacitanceSample> {
        let capacitance = problem.capacitance(
            &self.terminals[0],
            &self.terminals[1],
            self.controls,
            cancelled,
        )?;
        let field = capacitance.unit_voltage_field();
        // A displacement of one coordinate metre defines the derivative;
        // virtual_work is a differential, not a finite one-metre mesh update.
        let derivative_f_per_m = finite(2.0 * field.virtual_work(&self.motion, cancelled)?)?;
        let mut charge_defect_per_v = field
            .charge_balance_residual_c()
            .abs()
            .max(field.maximum_free_charge_residual_c());
        for (terminal, sign) in self.terminals.iter().zip([1.0, -1.0]) {
            let mut charge = Dd::ZERO;
            for &vertex in terminal {
                check(cancelled)?;
                let value = field.nodal_boundary_charge_c()[vertex]
                    .ok_or(Error::Invalid("original electrode reaction is absent"))?;
                charge = charge + Dd::from_f64(value);
            }
            charge_defect_per_v = charge_defect_per_v.max(
                finite((charge - Dd::from_f64(sign * capacitance.capacitance_f())).to_f64())?.abs(),
            );
        }
        let sample = CapacitanceSample {
            coordinate_m: coordinate,
            capacitance_f: capacitance.capacitance_f(),
            derivative_f_per_m,
            charge_defect_per_v,
            energy_defect_per_v2: field.energy_balance_residual_j().abs(),
        };
        check(cancelled)?;
        Ok(sample)
    }

    pub(super) fn admit_amplitude(&self, sample: CapacitanceSample, voltage: f64) -> Result<()> {
        let charge_defect = finite(voltage.abs() * sample.charge_defect_per_v)?;
        let energy_defect = finite(voltage * (voltage * sample.energy_defect_per_v2))?.abs();
        if charge_defect > self.controls.charge_tolerance_c
            || energy_defect > self.controls.energy_tolerance_j
        {
            Err(Error::Invalid(
                "actual-voltage field charge or energy residual exceeds its original allowance",
            ))
        } else {
            Ok(())
        }
    }
}
