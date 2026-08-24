// Photo mode (presentation-plane only): one body class toggled by the
// P key. The class hides every HUD surface and applies a period plate
// grade + frame border to the renderer canvas. I exports the same grade,
// deterministic grain, and border into PNG pixels; the physics/render
// pipelines remain untouched.

export interface PhotoModeToggle {
  /** True after the toggle if photo mode is now ON. */
  readonly active: boolean;
}

const PHOTO_CLASS = "wf-photo";

/** 3840×2160. The export path refuses larger readbacks before allocating. */
export const PHOTO_MAX_PIXELS = 8_294_400;

export const PHOTO_FILTER = "sepia(55%) contrast(106%) brightness(102%) saturate(85%)";

export type PhotoExportRefusalCode =
  | "photo-canvas-missing"
  | "photo-canvas-empty"
  | "photo-pixel-budget-exceeded"
  | "photo-identity-missing"
  | "photo-context-unavailable"
  | "photo-readback-failed"
  | "photo-encode-failed"
  | "photo-download-failed";

export interface PhotoExportRefusal {
  readonly code: PhotoExportRefusalCode;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type PhotoExportResult =
  | {
      readonly ok: true;
      readonly filename: string;
      readonly width: number;
      readonly height: number;
      readonly bytes: number;
    }
  | { readonly ok: false; readonly refusal: PhotoExportRefusal };

type PhotoExportFailure = Extract<PhotoExportResult, { readonly ok: false }>;

export interface PhotoUrlApi {
  createObjectURL(blob: Blob): string;
  revokeObjectURL(url: string): void;
}

function refused(
  code: PhotoExportRefusalCode,
  message: string,
  ...rankedRepairs: string[]
): PhotoExportFailure {
  return { ok: false, refusal: { code, message, rankedRepairs } };
}

/** Bounded size admission, exported so cap and cap+1 stay cheap to test. */
export function admitPhotoSize(
  width: number,
  height: number,
): { readonly ok: true; readonly pixels: number } | { readonly ok: false; readonly refusal: PhotoExportRefusal } {
  if (!Number.isSafeInteger(width) || !Number.isSafeInteger(height) || width <= 0 || height <= 0) {
    return refused(
      "photo-canvas-empty",
      `photo canvas has invalid dimensions ${width}×${height}`,
      "wait for the renderer to publish a non-empty frame",
    );
  }
  const pixels = width * height;
  if (!Number.isSafeInteger(pixels) || pixels > PHOTO_MAX_PIXELS) {
    return refused(
      "photo-pixel-budget-exceeded",
      `photo canvas requests ${pixels} pixels; limit is ${PHOTO_MAX_PIXELS}`,
      "reduce the browser viewport before exporting",
      "capture at 3840×2160 or smaller",
    );
  }
  return { ok: true, pixels };
}

/** Stable 32-bit seed from the run identity (FNV-1a). */
export function photoGrainSeed(identity: string): number {
  let hash = 0x811c9dc5;
  for (let i = 0; i < identity.length; i += 1) {
    hash ^= identity.charCodeAt(i);
    hash = Math.imul(hash, 0x01000193);
  }
  return hash >>> 0;
}

function grainDelta(seed: number, pixelIndex: number): number {
  let word = (seed ^ Math.imul(pixelIndex + 1, 0x9e3779b1)) >>> 0;
  word ^= word >>> 16;
  word = Math.imul(word, 0x7feb352d);
  word ^= word >>> 15;
  word = Math.imul(word, 0x846ca68b);
  word ^= word >>> 16;
  return ((word >>> 27) & 31) - 15;
}

/** Apply identity-keyed monochrome plate grain in place. Alpha is preserved. */
export function applyPhotoGrain(bytes: Uint8ClampedArray, seed: number): void {
  const pixels = Math.floor(bytes.length / 4);
  for (let pixel = 0; pixel < pixels; pixel += 1) {
    const offset = pixel * 4;
    const delta = grainDelta(seed, pixel);
    bytes[offset] = (bytes[offset] ?? 0) + delta;
    bytes[offset + 1] = (bytes[offset + 1] ?? 0) + delta;
    bytes[offset + 2] = (bytes[offset + 2] ?? 0) + delta;
  }
}

export function photoFilename(identity: string): string | null {
  const safe = identity
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9_-]+/g, "-")
    .replace(/^-+|-+$/g, "")
    .slice(0, 64);
  return safe.length > 0 ? `wright-flyer-${safe}.png` : null;
}

/**
 * Bake and download an identity-bound PNG from the live renderer canvas.
 * CSS alone is not accepted as export evidence: filter, grain, and border
 * are all applied to the offscreen pixels before encoding.
 */
export async function exportInstantPhoto(
  root: Document,
  source: HTMLCanvasElement | null,
  identity: string,
  urls: PhotoUrlApi = URL,
): Promise<PhotoExportResult> {
  if (source === null) {
    return refused(
      "photo-canvas-missing",
      "the renderer has not published a canvas",
      "start a flight or the attract-mode renderer before exporting",
    );
  }
  const admitted = admitPhotoSize(source.width, source.height);
  if (!admitted.ok) {
    return admitted;
  }
  const filename = photoFilename(identity);
  if (filename === null) {
    return refused(
      "photo-identity-missing",
      "the current run has no export identity",
      "wait for the simulation ready event before exporting",
    );
  }

  const plate = root.createElement("canvas");
  plate.width = source.width;
  plate.height = source.height;
  const ctx = plate.getContext("2d");
  if (ctx === null) {
    return refused(
      "photo-context-unavailable",
      "the browser could not create a 2-D plate canvas",
      "enable browser canvas support",
      "retry without a GPU-deny policy",
    );
  }

  try {
    ctx.filter = PHOTO_FILTER;
    ctx.drawImage(source, 0, 0, plate.width, plate.height);
    ctx.filter = "none";
    const image = ctx.getImageData(0, 0, plate.width, plate.height);
    applyPhotoGrain(image.data, photoGrainSeed(identity));
    ctx.putImageData(image, 0, 0);
    const shortest = Math.min(plate.width, plate.height);
    const inset = Math.min(
      Math.max(2, Math.round(shortest * 0.015)),
      Math.max(0, Math.floor((shortest - 1) / 2)),
    );
    ctx.strokeStyle = "rgba(243, 233, 210, 0.9)";
    ctx.lineWidth = Math.max(1, Math.round(inset * 0.18));
    ctx.strokeRect(inset, inset, plate.width - 2 * inset, plate.height - 2 * inset);
  } catch (error) {
    return refused(
      "photo-readback-failed",
      `the rendered frame could not be read: ${error instanceof Error ? error.message : String(error)}`,
      "keep image assets on the same origin",
      "retry after the renderer publishes a complete frame",
    );
  }

  let blob: Blob | null;
  try {
    blob = await new Promise<Blob | null>((resolve) => plate.toBlob(resolve, "image/png"));
  } catch (error) {
    return refused(
      "photo-encode-failed",
      `PNG encoding failed: ${error instanceof Error ? error.message : String(error)}`,
      "retry in a browser with PNG canvas encoding",
    );
  }
  if (blob === null) {
    return refused(
      "photo-encode-failed",
      "the browser returned no PNG bytes",
      "retry in a browser with PNG canvas encoding",
    );
  }

  let objectUrl: string | null = null;
  try {
    objectUrl = urls.createObjectURL(blob);
    const anchor = root.createElement("a");
    anchor.href = objectUrl;
    anchor.download = filename;
    anchor.click();
  } catch (error) {
    return refused(
      "photo-download-failed",
      `the PNG was encoded but download failed: ${error instanceof Error ? error.message : String(error)}`,
      "allow downloads from this page",
      "retry from a top-level browser tab",
    );
  } finally {
    if (objectUrl !== null) {
      try {
        urls.revokeObjectURL(objectUrl);
      } catch {
        // Download already fired; failure to revoke must not turn a
        // delivered PNG into a false failure or an unhandled rejection.
      }
    }
  }

  return { ok: true, filename, width: plate.width, height: plate.height, bytes: blob.size };
}

/** Toggle photo mode on a document root. Returns the new state. Pure
 * DOM-class mutation; safe to call repeatedly. */
export function togglePhotoMode(root: Document): PhotoModeToggle {
  const active = !root.body.classList.contains(PHOTO_CLASS);
  root.body.classList.toggle(PHOTO_CLASS, active);
  return { active };
}

/** Current state without toggling. */
export function photoModeActive(root: Document): boolean {
  return root.body.classList.contains(PHOTO_CLASS);
}
