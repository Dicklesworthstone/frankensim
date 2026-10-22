//! Spatial pressure traction on the ORIGINAL material blade eigenbasis.
//!
//! Each port supplies p(x,t) = psi(x) P(t), where psi is an explicit signed,
//! dimensionless P1 trace and P is in pascals. Its generalized force is
//! f_i = P integral(phi_i psi dA); its conjugate volume flow is the SAME
//! overlap dotted with mean modal velocity. No net-force normalization is
//! performed: a zero-resultant pressure field can still bend a blade.
//!
//! These are additional prescribed external loads, not a solved channel
//! pressure field. They do not replace the uniform bore pressure, change its
//! swept-volume projection, or add another flow to either acoustic inlet.

use super::{AcousticRealizeError, ApertureDrive, ModalAperture, ModalApertureFrame, invalid};
use super::material::MaterialModalPlate;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;

/// Explicit spatial pressure trace bound to one original material specimen.
/// Values are dimensionless and unnormalized; an amplitude P applies psi*P Pa.
/// The source label identifies the prescribed model or independent measurement.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialPressureShape {
    specimen: ContentHash,
    nodal_shape: Vec<f64>,
    provenance: String,
}
impl MaterialPressureShape {
    /// Bind values in original node order. The receiving valve subsequently
    /// checks source identity, complete dimensions, finiteness and work budget.
    pub fn from_plate(plate: &MaterialModalPlate, nodal_shape: Vec<f64>, provenance: String) -> Self {
        Self { specimen: plate.source().specimen_identity(), nodal_shape, provenance }
    }
    /// Exact geometry/material source identity; the consumer never substitutes it.
    pub const fn specimen_identity(&self) -> ContentHash { self.specimen }
    /// Original dimensionless P1 values. No averaging, clipping or rescaling.
    pub fn nodal_shape(&self) -> &[f64] { &self.nodal_shape }
    /// Independent source or explicit authored-model label.
    pub fn provenance(&self) -> &str { &self.provenance }
}

/// Pressure and its conjugate projected volume flow over one accepted sample.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaterialPressureFrame {
    /// Held reference amplitude [Pa]; local pressure is this times psi(x).
    pub amplitude_pa: f64,
    /// integral(psi * mean blade velocity dA) [m³/s], positive toward opening.
    /// This generalized pressure-port flow is NOT an extra duct through-flow.
    pub flow_m3_s: f64,
    /// Signed external work P*Q*dt [J], already included in mechanical work.
    pub work_j: f64,
}

/// Retained original trace and its pressure/velocity overlaps [m²].
#[derive(Clone, Debug)]
pub struct MaterialPressurePort {
    shape: MaterialPressureShape,
    modal_area_m2: Vec<f64>,
}
impl MaterialPressurePort {
    /// Immutable original pressure data and specimen identity.
    pub const fn shape(&self) -> &MaterialPressureShape { &self.shape }
    /// One signed overlap per complete valve coordinate; other blades are zero.
    pub fn modal_area_m2(&self) -> &[f64] { &self.modal_area_m2 }
    pub(super) fn observe(&self, amplitude_pa: f64, velocity: &[f64], dt: f64) -> MaterialPressureFrame {
        let flow_m3_s = self.modal_area_m2.iter().zip(velocity).map(|(a, v)| a*v).sum::<f64>();
        MaterialPressureFrame { amplitude_pa, flow_m3_s, work_j: dt*amplitude_pa*flow_m3_s }
    }
}

impl ModalAperture {
    /// Attach one additional pressure trace to an actual material blade before
    /// stepping. Consistent P1 triangle quadrature comes from the shared plate
    /// pressure owner and projects the already retained displacement columns.
    /// Elastic modes, full viscosity, Maxwell histories and contact stay intact.
    /// Signed/zero-resultant fields are allowed; an identically ineffective port
    /// refuses instead of claiming that it can drive these retained modes.
    ///
    /// This is an independently prescribed load. Reusing the already coupled
    /// uniform bore-pressure difference here would count the same source twice.
    pub fn with_material_pressure_shape(mut self, blade: usize, shape: MaterialPressureShape,
        max_projection_terms: usize, gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        if self.samples != 0 || blade >= 2 || self.pressure_ports[blade].is_some() {
            return Err(invalid("material pressure requires one cold attachment per original blade"));
        }
        let plate = self.material_plates.get(blade)
            .ok_or_else(|| invalid("material pressure names an absent original blade"))?;
        if shape.specimen != plate.source().specimen_identity()
            || shape.provenance.trim().is_empty() || shape.provenance.len() > 4096
            || shape.nodal_shape.len() != plate.displacement_shapes().len()
            || shape.nodal_shape.iter().any(|p| !p.is_finite()) {
            return Err(invalid("material pressure requires the same specimen, complete finite nodal trace and source"));
        }
        let terms = shape.nodal_shape.len().checked_mul(plate.mode_range().len())
            .ok_or_else(|| invalid("material pressure projection work overflow"))?;
        let remaining = max_projection_terms.checked_sub(terms)
            .ok_or_else(|| invalid("material pressure projection exceeds its work budget"))?;
        let loads = plate.source().nodal_pressure_loads(&[shape.nodal_shape.clone()], remaining)?;
        let mut area = vec![0.0; self.state.len()];
        for (local, mode) in plate.mode_range().enumerate() {
            area[mode] = plate.displacement_shapes().iter().zip(&loads[0])
                .map(|(row, load)| row[local]*load).sum();
        }
        if area.iter().any(|a| !a.is_finite()) || area.iter().all(|a| *a == 0.0) {
            return Err(invalid("material pressure has no finite nonzero retained overlap"));
        }
        if gate.is_requested() { return Err(AcousticRealizeError::Cancelled); }
        if self.pressure_ports.iter().all(Option::is_none) {
            for name in ["pressure/0", "pressure/1"] {
                self.ports.schemas.push(crate::port_adoption::schema(
                    "modal-aperture", name, fs_ports::PortKind::FluidPressureFlux));
                self.ports.charts.push(fs_ports::PortCoordinateChart::CanonicalSi);
            }
        }
        self.pressure_ports[blade] = Some(MaterialPressurePort { shape, modal_area_m2: area });
        Ok(self)
    }

    /// Bound independent pressure ports in original blade order.
    pub const fn material_pressure_ports(&self) -> &[Option<MaterialPressurePort>; 2] {
        &self.pressure_ports
    }

    /// Refuse nonfinite loads or nonzero pressure addressed to an absent port.
    pub fn validate_pressures(&self, pressures_pa: [f64; 2]) -> Result<(), AcousticRealizeError> {
        if pressures_pa.iter().enumerate().any(|(i, p)| !p.is_finite()
            || (*p != 0.0 && self.pressure_ports[i].is_none())) {
            return Err(invalid("spatial pressure needs finite pascals and an existing material port"));
        }
        Ok(())
    }

    /// Advance spatial traction, uniform acoustic pressure, contact and all
    /// material states through the same existing implicit structural solve.
    /// A refusal leaves every modal and hereditary coordinate unchanged.
    pub fn step_with_pressures(&mut self, drive: ApertureDrive, forces_n: [f64; 2],
        pressures_pa: [f64; 2]) -> Result<ModalApertureFrame, AcousticRealizeError>
    {
        let frame = self.preview_with_pressures(drive, forces_n, pressures_pa)?;
        self.accept(&frame);
        Ok(frame)
    }
}
