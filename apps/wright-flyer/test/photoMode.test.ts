// Instant-photo battery (E9.2a): bounded admission, deterministic
// identity-keyed grain, baked grade/border, typed encoder failures,
// and the final browser download seam.
// Repro: node --test test/photoMode.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  PHOTO_FILTER,
  PHOTO_MAX_PIXELS,
  admitPhotoSize,
  applyPhotoGrain,
  exportInstantPhoto,
  photoFilename,
  photoGrainSeed,
  type PhotoUrlApi,
} from "../src/photoMode.ts";

function fakeSource(width: number, height: number): HTMLCanvasElement {
  return { width, height } as HTMLCanvasElement;
}

interface PlateHarness {
  readonly plate: HTMLCanvasElement;
  readonly filtersAtDraw: string[];
  readonly strokes: number[];
  readonly putBytes: Uint8ClampedArray[];
}

function fakePlate(
  width: number,
  height: number,
  encoded: Blob | null = new Blob(["png"], { type: "image/png" }),
): PlateHarness {
  const filtersAtDraw: string[] = [];
  const strokes: number[] = [];
  const putBytes: Uint8ClampedArray[] = [];
  const initial = new Uint8ClampedArray(width * height * 4).fill(100);
  for (let i = 3; i < initial.length; i += 4) {
    initial[i] = 255;
  }
  const context = {
    filter: "none",
    strokeStyle: "",
    lineWidth: 0,
    drawImage(): void {
      filtersAtDraw.push(this.filter);
    },
    getImageData(): ImageData {
      return { data: initial.slice(), width, height, colorSpace: "srgb" } as ImageData;
    },
    putImageData(image: ImageData): void {
      putBytes.push(image.data.slice());
    },
    strokeRect(): void {
      strokes.push(this.lineWidth);
    },
  };
  const plate = {
    width: 0,
    height: 0,
    getContext(kind: string): CanvasRenderingContext2D | null {
      return kind === "2d" ? (context as unknown as CanvasRenderingContext2D) : null;
    },
    toBlob(callback: BlobCallback): void {
      callback(encoded);
    },
  } as unknown as HTMLCanvasElement;
  return { plate, filtersAtDraw, strokes, putBytes };
}

function fakeDocument(
  plate: HTMLCanvasElement,
  downloads: Array<{ href: string; filename: string }>,
): Document {
  return {
    createElement(tag: string): HTMLCanvasElement | HTMLAnchorElement {
      if (tag === "canvas") {
        return plate;
      }
      if (tag === "a") {
        const anchor = {
          href: "",
          download: "",
          click(): void {
            downloads.push({ href: this.href, filename: this.download });
          },
        };
        return anchor as HTMLAnchorElement;
      }
      throw new Error(`unexpected element ${tag}`);
    },
  } as unknown as Document;
}

function fakeUrls(events: string[]): PhotoUrlApi {
  return {
    createObjectURL(blob): string {
      events.push(`create:${blob.size}`);
      return "blob:wf-photo";
    },
    revokeObjectURL(url): void {
      events.push(`revoke:${url}`);
    },
  };
}

test("size admission accepts the 4K cap and refuses cap+1 before allocation", () => {
  assert.deepEqual(admitPhotoSize(3840, 2160), { ok: true, pixels: PHOTO_MAX_PIXELS });
  const over = admitPhotoSize(PHOTO_MAX_PIXELS + 1, 1);
  assert.equal(over.ok, false);
  if (!over.ok) {
    assert.equal(over.refusal.code, "photo-pixel-budget-exceeded");
  }
  for (const [width, height] of [[0, 10], [10, 0], [1.5, 2], [Number.NaN, 2]]) {
    const empty = admitPhotoSize(width, height);
    assert.equal(empty.ok, false);
    if (!empty.ok) {
      assert.equal(empty.refusal.code, "photo-canvas-empty");
    }
  }
});

test("identity-keyed plate grain is deterministic, non-vacuous, and alpha-safe", () => {
  const seed = photoGrainSeed("run-1903");
  assert.equal(seed, photoGrainSeed("run-1903"));
  assert.notEqual(seed, photoGrainSeed("run-1904"));
  const a = new Uint8ClampedArray(64).fill(120);
  const b = a.slice();
  for (let i = 3; i < a.length; i += 4) {
    a[i] = 255;
    b[i] = 255;
  }
  const original = a.slice();
  applyPhotoGrain(a, seed);
  applyPhotoGrain(b, seed);
  assert.deepEqual(a, b);
  assert.notDeepEqual(a, original, "grain changes RGB pixels");
  for (let i = 3; i < a.length; i += 4) {
    assert.equal(a[i], 255, "grain must preserve alpha");
  }
});

test("successful export bakes grade, grain, border, and downloads an identity filename", async () => {
  const harness = fakePlate(2, 2);
  const downloads: Array<{ href: string; filename: string }> = [];
  const urlEvents: string[] = [];
  const result = await exportInstantPhoto(
    fakeDocument(harness.plate, downloads),
    fakeSource(2, 2),
    "Run Intent: 1903/A",
    fakeUrls(urlEvents),
  );
  assert.deepEqual(result, {
    ok: true,
    filename: "wright-flyer-run-intent-1903-a.png",
    width: 2,
    height: 2,
    bytes: 3,
  });
  assert.deepEqual(harness.filtersAtDraw, [PHOTO_FILTER]);
  assert.equal(harness.putBytes.length, 1, "grain pixels are written back");
  assert.equal(harness.strokes.length, 1, "plate border is baked");
  assert.deepEqual(downloads, [
    { href: "blob:wf-photo", filename: "wright-flyer-run-intent-1903-a.png" },
  ]);
  assert.deepEqual(urlEvents, ["create:3", "revoke:blob:wf-photo"]);
});

test("missing canvas, missing identity, and encoder null fail closed without download", async () => {
  const downloads: Array<{ href: string; filename: string }> = [];
  const urls: string[] = [];
  const harness = fakePlate(2, 2, null);
  const root = fakeDocument(harness.plate, downloads);

  const missing = await exportInstantPhoto(root, null, "run", fakeUrls(urls));
  assert.equal(missing.ok, false);
  if (!missing.ok) assert.equal(missing.refusal.code, "photo-canvas-missing");

  const unbound = await exportInstantPhoto(root, fakeSource(2, 2), "  ", fakeUrls(urls));
  assert.equal(unbound.ok, false);
  if (!unbound.ok) assert.equal(unbound.refusal.code, "photo-identity-missing");

  const failed = await exportInstantPhoto(root, fakeSource(2, 2), "run", fakeUrls(urls));
  assert.equal(failed.ok, false);
  if (!failed.ok) assert.equal(failed.refusal.code, "photo-encode-failed");
  assert.deepEqual(downloads, []);
  assert.deepEqual(urls, []);
});

test("filename normalization is stable and bounded", () => {
  assert.equal(photoFilename(" Run Intent: 1903/A "), "wright-flyer-run-intent-1903-a.png");
  assert.equal(photoFilename("   "), null);
  assert.ok(photoFilename("x".repeat(100))!.length <= "wright-flyer-.png".length + 64);
});
