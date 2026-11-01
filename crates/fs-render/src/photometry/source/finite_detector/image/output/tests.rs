use super::super::super::{FiniteDetectorSpectrum, tests::with_cx};
use super::super::DetectorPixel;
use super::*;
use crate::instances::RigidTransform;
use crate::photometry::observer::xyz_at_bin;
use crate::photometry::source::{BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::measurement::raster::DetectorRaster;

const BIN: usize = 31;
const SAMPLES: u32 = 30;

fn stokes(intensity: f64) -> Stokes {
    Stokes {
        i: intensity,
        q: -0.5 * intensity,
        u: 0.25 * intensity,
        v: 0.0,
    }
}

fn pixel(watts: f64) -> DetectorPixel {
    // A monochromatic 15 W packet is detected on every one of 30 launches.
    // Each cell receives 2*watts of those launches, and the fifth cell is the
    // explicit crop. This gives a known nontrivial mean and Bernoulli error.
    let fraction = watts / 15.0;
    let error = 15.0 * (fraction * (1.0 - fraction) / f64::from(SAMPLES - 1)).sqrt();
    let xyz = xyz_at_bin(BIN).map(|value| LUMINOUS_EFFICACY_LM_W * value);
    DetectorPixel {
        radiant: stokes(watts),
        photopic: stokes(watts * xyz[1]),
        xz: [watts * xyz[0], watts * xyz[2]],
        radiant_error: [error, 0.5 * error, 0.25 * error, 0.0],
        photopic_error: [
            error * xyz[1],
            0.5 * error * xyz[1],
            0.25 * error * xyz[1],
            0.0,
        ],
        xz_error: [error * xyz[0], error * xyz[2]],
        received_samples: (2.0 * watts) as u32,
        detected_paths: (2.0 * watts) as u64,
        max_contact_resolution_m: 1e-7,
    }
}

fn image() -> DetectorImage {
    let xyz = xyz_at_bin(BIN).map(|value| LUMINOUS_EFFICACY_LM_W * value);
    let mut spectral = [Stokes::unpolarized(0.0); SPECTRAL_BINS];
    spectral[BIN] = stokes(15.0 / BIN_WIDTH_NM);
    DetectorImage {
        raster: DetectorRaster::try_new(RigidTransform::identity(), 2.0, 1.0, 2, 2).unwrap(),
        pixels: vec![pixel(1.0), pixel(2.0), pixel(3.0), pixel(4.0)],
        cropped: pixel(5.0),
        total: FiniteDetectorSpectrum {
            spectral,
            spectral_error: [[0.0; 4]; SPECTRAL_BINS],
            radiant: stokes(15.0),
            photopic: stokes(15.0 * xyz[1]),
            radiant_error: [0.0; 4],
            photopic_error: [0.0; 4],
            samples: SAMPLES,
            detected_paths: u64::from(SAMPLES),
            received_samples: SAMPLES,
            max_contact_resolution_m: 1e-7,
        },
    }
}

fn channel<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a [f32] {
    &decoded
        .channels
        .iter()
        .find(|channel| channel.name == name)
        .unwrap()
        .data
}

fn attribute<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a ExrAttribute {
    decoded
        .attributes
        .iter()
        .find(|attribute| attribute.name == name)
        .unwrap()
}

fn double_attribute(decoded: &fs_img::DecodedExr, name: &str) -> f64 {
    let value = attribute(decoded, name);
    assert_eq!(value.ty, "double");
    f64::from_le_bytes(value.value.as_slice().try_into().unwrap())
}

#[test]
fn g1_openexr_roundtrip_retains_absolute_pixels_signed_stokes_and_crop_metadata() {
    with_cx(|cx, _| {
        let image = image();
        let bytes = image.to_openexr(cx, ExrWriteLimits::UNBOUNDED).unwrap();
        assert_eq!(
            bytes,
            image.to_openexr(cx, ExrWriteLimits::UNBOUNDED).unwrap()
        );
        let decoded = fs_img::read_exr(&bytes).unwrap();
        assert_eq!((decoded.width, decoded.height), (2, 2));
        assert_eq!(decoded.channels.len(), 20);
        assert!(
            decoded
                .channels
                .iter()
                .all(|channel| channel.ty == PixelType::Float)
        );
        // The 0.5 m² pixels retain top-to-bottom row order and absolute scale.
        assert_eq!(channel(&decoded, "irradiance.I"), &[2.0, 4.0, 6.0, 8.0]);
        assert_eq!(channel(&decoded, "irradiance.Q"), &[-1.0, -2.0, -3.0, -4.0]);
        assert_eq!(channel(&decoded, "irradiance.U"), &[0.5, 1.0, 1.5, 2.0]);
        for (index, pixel) in image.pixels().iter().enumerate() {
            assert_eq!(
                channel(&decoded, "illuminance.I")[index],
                (2.0 * pixel.photopic_stokes_lm().i) as f32
            );
            assert_eq!(
                channel(&decoded, "cie.X")[index],
                (2.0 * pixel.xyz_flux()[0]) as f32
            );
            assert_eq!(
                channel(&decoded, "cie.Z")[index],
                (2.0 * pixel.xyz_flux()[2]) as f32
            );
            assert_eq!(
                channel(&decoded, "stderr.irradiance.I")[index],
                (2.0 * pixel.radiant_standard_error_w()[0]) as f32
            );
            assert_eq!(
                channel(&decoded, "stderr.illuminance.Q")[index],
                (2.0 * pixel.photopic_standard_error_lm()[1]) as f32
            );
            assert_eq!(
                channel(&decoded, "stderr.cie.X")[index],
                (2.0 * pixel.xyz_standard_error()[0]) as f32
            );
        }
        assert_eq!(double_attribute(&decoded, "frankensim.pixelArea_m2"), 0.5);
        assert_eq!(double_attribute(&decoded, "frankensim.pixelWidth_m"), 1.0);
        assert_eq!(double_attribute(&decoded, "frankensim.pixelHeight_m"), 0.5);
        assert_eq!(double_attribute(&decoded, "frankensim.cropped_W"), 5.0);
        assert_eq!(double_attribute(&decoded, "frankensim.total_W"), 15.0);
        assert_eq!(
            double_attribute(&decoded, "frankensim.contactResolution_m"),
            1e-7
        );
        assert_eq!(attribute(&decoded, "frankensim.sensorSamples").value, b"30");
        let units =
            std::str::from_utf8(&attribute(&decoded, "frankensim.sensorUnits").value).unwrap();
        assert!(
            units.contains("W/m^2")
                && units.contains("lux")
                && units.contains("cie.Y=illuminance.I")
        );
        let rows =
            std::str::from_utf8(&attribute(&decoded, "frankensim.sensorRows").value).unwrap();
        assert!(rows.contains("row0=top"));
        let spectrum =
            std::str::from_utf8(&attribute(&decoded, "frankensim.sensorSpectrum").value).unwrap();
        assert!(spectrum.contains("80 midpoint bins") && spectrum.contains("CIE 1931"));
    });
}

#[test]
fn g4_openexr_checks_exact_byte_limits_before_pixel_projection_and_obeys_cancellation() {
    with_cx(|cx, gate| {
        let mut image = image();
        let attributes = attributes(&image);
        let layout = CHANNELS.map(|name| (name, PixelType::Float));
        let requirements =
            fs_img::exr_write_requirements_for_layout(2, 2, &layout, &attributes).unwrap();
        let exact = ExrWriteLimits {
            max_scratch_bytes: requirements.scratch_bytes,
            max_output_bytes: requirements.output_bytes,
        };
        assert_eq!(
            image.to_openexr(cx, exact).unwrap().len() as u64,
            requirements.output_bytes
        );
        assert!(matches!(
            image.to_openexr(
                cx,
                ExrWriteLimits {
                    max_scratch_bytes: exact.max_scratch_bytes - 1,
                    ..exact
                }
            ),
            Err(SensorImageExportError::Image(ImgError::ResourceLimit {
                resource: "EXR ordering scratch",
                ..
            }))
        ));
        // This otherwise valid binary64 measurement cannot fit FLOAT. The
        // output ceiling still refuses first, before pixel-plane allocation.
        image.pixels[0].radiant = Stokes::unpolarized(1e100);
        assert!(matches!(
            image.to_openexr(
                cx,
                ExrWriteLimits {
                    max_output_bytes: exact.max_output_bytes - 1,
                    ..exact
                }
            ),
            Err(SensorImageExportError::Image(ImgError::ResourceLimit {
                resource: "EXR encoded output",
                ..
            }))
        ));
        gate.request();
        assert_eq!(
            image.to_openexr(cx, exact),
            Err(SensorImageExportError::Cancelled)
        );
    });
}

#[test]
fn g4_openexr_refuses_float_overflow_and_nonzero_underflow_without_clamping() {
    with_cx(|cx, _| {
        let mut image = image();
        for power in [1e100, 1e-200] {
            image.pixels[0].radiant = Stokes::unpolarized(power);
            assert_eq!(
                image.to_openexr(cx, ExrWriteLimits::UNBOUNDED),
                Err(SensorImageExportError::UnrepresentableSample {
                    pixel: 0,
                    channel: "irradiance.I"
                }),
            );
        }
    });
}
