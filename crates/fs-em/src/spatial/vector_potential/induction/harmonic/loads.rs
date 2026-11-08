//! Consistent mean and twice-frequency Lorentz loads from one native AC field.
//!
//! For peak phasors with `exp(-i omega t)`, the physical force is
//! `f(t) = f_mean + Re(f_2 exp(-2 i omega t))`, where
//! `f_mean = Re(J cross conj(B))/2` and `f_2 = (J cross B)/2`.
//! Keeping only the mean loses the force that drives magnetic vibration.
//! Here J is P1 and B is constant in each original tetrahedron, so all three
//! real coefficient fields are P1. The original consistent load integrator
//! preserves their force, torque and virtual work without vertex lumping.
//!
//! These are conduction-current volume forces. The static-bias child retains
//! fundamental forcing in stationary, DC-current-free selected conductors.
//! General moving/DC-driven cross terms, impressed winding-current forces,
//! magnetization, magnetostriction and interface traction are outside this transfer. Do not add an overlapping enclosing
//! Maxwell-stress load to the same body. Mechanical feedback, geometry changes
//! and vibration-induced electrical sidebands require a coupled evolution;
//! these fixed-geometry coefficients supply one-way mechanical forcing only.

use super::super::mechanical::{
    SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads,
    load::{CellLoad, integrate},
};
use super::{SpatialHarmonicInductionResponse, cross};
use fs_material::conductor::spatial::SpatialConductorMaterial;

mod bias;
pub use bias::SpatialHarmonicBiasLorentzLoads;
mod cross;
pub use cross::SpatialHarmonicLorentzCrossLoads;

/// Mechanical Fourier coefficients in the original field coordinates.
///
/// Every coefficient uses original mesh vertex order and the existing load
/// owner's receiving-frame wrench and P1 virtual-power APIs. In particular,
/// the real and imaginary second-harmonic parts are PEAK coefficients, not
/// RMS values or additional cycle means. At electromagnetic phase `theta`,
/// nodal forces and wrenches reconstruct as `mean + real*cos(2 theta) +
/// imaginary*sin(2 theta)` for the native negative-exponential convention.
#[derive(Debug)]
pub struct SpatialHarmonicLorentzLoads {
    second_harmonic_omega_rad_s: f64,
    mean: SpatialLorentzLoads,
    second_harmonic_real: SpatialLorentzLoads,
    second_harmonic_imaginary: SpatialLorentzLoads,
}

impl SpatialHarmonicLorentzLoads {
    /// Angular frequency [rad/s] of the oscillating force, twice the EM carrier.
    pub fn second_harmonic_omega_rad_s(&self) -> f64 {
        self.second_harmonic_omega_rad_s
    }

    /// Cycle-mean nodal forces, wrench and power against a real P1 velocity.
    ///
    /// Pairing with the field solve's prescribed steady Eulerian velocity
    /// recovers the selected cells' mean Lorentz-mechanical power. Pairing with
    /// another velocity is virtual work, not a closed field/mechanical solve.
    pub fn mean(&self) -> &SpatialLorentzLoads {
        &self.mean
    }

    /// Real part [N, N m] of the PEAK complex force/wrench at twice frequency.
    pub fn second_harmonic_real(&self) -> &SpatialLorentzLoads {
        &self.second_harmonic_real
    }

    /// Imaginary part [N, N m] of the PEAK complex force/wrench at twice frequency.
    /// It multiplies positive `sin(2 omega t)` in real time reconstruction.
    pub fn second_harmonic_imaginary(&self) -> &SpatialLorentzLoads {
        &self.second_harmonic_imaginary
    }
}

impl<C: SpatialConductorMaterial> SpatialHarmonicInductionResponse<'_, '_, C> {
    /// Extract mean AND twice-frequency conduction-current loads on a body.
    ///
    /// Selection names distinct original conducting tetrahedra in any order.
    /// It is canonicalized before summation; unselected vertices receive zero
    /// force. The geometry, J phasors and B phasors come from this response's
    /// original source-bound field, without a static reset or another solve.
    /// Each coefficient reuses exact tetrahedral P1 mass weights, compensated
    /// scatter and the same local-reference torque reduction as transient
    /// induction. Integration costs O(mesh vertices + selected cells), with
    /// three retained real nodal-force arrays and bounded selected-cell data.
    /// `limits.max_vertices` bounds the original mesh, for each coefficient.
    ///
    /// Use `mean().wrench(...)` and the two second-harmonic `wrench(...)`
    /// methods for native rigid-body forcing, or their `nodal_forces_n()` for
    /// deformable mechanical receivers. All three must use the same placement
    /// and moment origin. Cancellation or refusal returns no partial loads
    /// and does not modify the field history.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicLorentzLoads, SpatialLorentzLoadError> {
        use SpatialLorentzLoadError::{Cancelled, Invalid};
        if cancelled() {
            return Err(Cancelled);
        }
        let problem = self.problem();
        let second_harmonic_omega_rad_s = 2.0 * self.omega_rad_s();
        if selected_cells.is_empty()
            || selected_cells.len() > limits.max_cells
            || problem.positions_m().len() > limits.max_vertices
            || !second_harmonic_omega_rad_s.is_finite()
        {
            return Err(Invalid(
                "empty, over-budget or unrepresentable harmonic Lorentz selection",
            ));
        }
        let mut selected = Vec::with_capacity(selected_cells.len());
        for &cell in selected_cells {
            if cancelled() {
                return Err(Cancelled);
            }
            let record = self
                .cells()
                .binary_search_by_key(&cell, |record| record.cell)
                .ok()
                .and_then(|index| self.cells().get(index))
                .ok_or(Invalid(
                    "harmonic Lorentz selection requires original conducting cells",
                ))?;
            let vertices = *problem.mesh().tets.get(cell).ok_or(Invalid(
                "harmonic Lorentz cell is outside its original mesh",
            ))?;
            let b = self.flux_density_t()[cell];
            let br = b.map(|v| v.re);
            let bi = b.map(|v| v.im);
            let mut real = [[0.0; 3]; 4];
            let mut imaginary = [[0.0; 3]; 4];
            for vertex in 0..4 {
                let jr = record.current_density_a_m2[vertex].map(|v| v.re);
                let ji = record.current_density_a_m2[vertex].map(|v| v.im);
                let rr = cross(jr, br);
                let ii = cross(ji, bi);
                let ri = cross(jr, bi);
                let ir = cross(ji, br);
                for axis in 0..3 {
                    real[vertex][axis] = f64::midpoint(rr[axis], -ii[axis]);
                    imaginary[vertex][axis] = f64::midpoint(ri[axis], ir[axis]);
                }
            }
            selected.push((cell, vertices, record, real, imaginary));
        }
        selected.sort_unstable_by_key(|entry| entry.0);
        if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(Invalid(
                "harmonic Lorentz selection contains a repeated cell",
            ));
        }
        let mean = integrate(
            problem.positions_m(),
            selected.iter().map(|(_, vertices, record, _, _)| CellLoad {
                vertices: *vertices,
                volume_m3: record.volume_m3,
                force_density_n_m3: &record.mean_lorentz_force_density_n_m3,
            }),
            limits,
            cancelled,
        )?;
        let second_harmonic_real = integrate(
            problem.positions_m(),
            selected
                .iter()
                .map(|(_, vertices, record, real, _)| CellLoad {
                    vertices: *vertices,
                    volume_m3: record.volume_m3,
                    force_density_n_m3: real,
                }),
            limits,
            cancelled,
        )?;
        let second_harmonic_imaginary = integrate(
            problem.positions_m(),
            selected
                .iter()
                .map(|(_, vertices, record, _, imaginary)| CellLoad {
                    vertices: *vertices,
                    volume_m3: record.volume_m3,
                    force_density_n_m3: imaginary,
                }),
            limits,
            cancelled,
        )?;
        Ok(SpatialHarmonicLorentzLoads {
            second_harmonic_omega_rad_s,
            mean,
            second_harmonic_real,
            second_harmonic_imaginary,
        })
    }
}
