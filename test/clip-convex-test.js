import assert from "node:assert";
import {it} from "vitest";
import {geoArea, geoAzimuthalEqualArea, geoClipConvex, geoPath, geoStream} from "d3-geo";
import {assertInDelta} from "./asserts.js";

const radians = Math.PI / 180;
const degrees = 180 / Math.PI;

function rad(ring) {
  return ring.map(([x, y]) => [x * radians, y * radians]);
}

function polygon(ring) {
  return {type: "Polygon", coordinates: [ring.concat([ring[0]])]};
}

// Streams `object` (spherical coordinates in radians) through the clip and
// collects the emitted polygons and lines (still in radians).
function clipCollect(object, clip) {
  const polygons = [], lines = [];
  let current = null, ring = null;
  geoStream(object, clip({
    polygonStart() { current = []; },
    polygonEnd() { polygons.push(current); current = null; },
    lineStart() { ring = []; },
    lineEnd() { if (current) current.push(ring); else lines.push(ring); ring = null; },
    point(x, y) { ring.push([x, y]); },
    sphere() {}
  }));
  return {polygons, lines};
}

function areaOf(rings) {
  return geoArea({type: "Polygon", coordinates: rings.map((ring) => {
    const closed = ring.map(([x, y]) => [x * degrees, y * degrees]);
    closed.push(closed[0].slice());
    return closed;
  })});
}

// Signed solid angle of the (small) spherical triangle with the given
// vertices (in radians), after Van Oosterom & Strackee (1983).
function sphericalTriangleArea(triangle) {
  const c = triangle.map(([lambda, phi]) => [
    Math.cos(phi) * Math.cos(lambda),
    Math.cos(phi) * Math.sin(lambda),
    Math.sin(phi)
  ]);
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return 2 * Math.atan2(
    dot(c[0], cross(c[1], c[2])),
    1 + dot(c[0], c[1]) + dot(c[1], c[2]) + dot(c[2], c[0])
  );
}

const worldRectangle = polygon(rad([[-180, -80], [180, -80], [180, 80], [-180, 80]]));
const triangle = [[-120, -60], [120, -60], [0, 60]];

it("clipConvex: clips the world rectangle to a triangle", () => {
  const clip = geoClipConvex(rad(triangle));
  const {polygons} = clipCollect(worldRectangle, clip);
  assert.strictEqual(polygons.length, 1);
  assert.strictEqual(polygons[0].length, 1);
  const expected = geoArea(polygon(triangle)); // the triangle is strictly inside the rectangle
  assertInDelta(areaOf(polygons[0]), expected, Math.abs(expected) * 1e-6);
  // Independently: the CCW ring covers everything but the small triangle.
  assertInDelta(areaOf(polygons[0]), 4 * Math.PI - Math.abs(sphericalTriangleArea(rad(triangle))), 1e-6);
});

it("clipConvex: preserves the input winding", () => {
  const reversed = [triangle[0], triangle[2], triangle[1]];
  const clip = geoClipConvex(rad(reversed));
  const {polygons} = clipCollect(worldRectangle, clip);
  assert.strictEqual(polygons.length, 1);
  const expected = geoArea(polygon(reversed));
  assertInDelta(areaOf(polygons[0]), expected, Math.abs(expected) * 1e-6);
  assertInDelta(areaOf(polygons[0]), Math.abs(sphericalTriangleArea(rad(reversed))), 1e-6);
});

it("clipConvex: a polygon fully inside is passed through identically", () => {
  const clip = geoClipConvex(rad(triangle));
  const inside = rad([[-10, -10], [10, -10], [0, 10]]);
  const {polygons} = clipCollect(polygon(inside), clip);
  assert.strictEqual(polygons.length, 1);
  assert.strictEqual(polygons[0].length, 1);
  assert.deepStrictEqual(polygons[0][0], inside);
});

it("clipConvex: the clip polygon itself is passed through identically", () => {
  const clip = geoClipConvex(rad([[-10, -10], [10, -10], [0, 10]]));
  const self = rad([[-10, -10], [10, -10], [0, 10]]);
  const {polygons} = clipCollect(polygon(self), clip);
  assert.strictEqual(polygons.length, 1);
  assert.deepStrictEqual(polygons[0][0], self);
});

it("clipConvex: a polygon fully outside produces an empty polygon", () => {
  const clip = geoClipConvex(rad(triangle));
  const outside = rad([[-10, -80], [10, -80], [0, -70]]);
  const {polygons} = clipCollect(polygon(outside), clip);
  assert.strictEqual(polygons.length, 0);
});

it("clipConvex: handles a clip polygon crossing the antimeridian", () => {
  // A viewport across the antimeridian, given with continuous longitudes.
  const viewport = rad([[170, -30], [190, -30], [190, -60], [170, -60]]);
  const clip = geoClipConvex(viewport);
  const {polygons} = clipCollect(worldRectangle, clip);
  assert.strictEqual(polygons.length, 1);
  const expected = geoArea(polygon([[170, -30], [-170, -30], [-170, -60], [170, -60]]));
  assertInDelta(areaOf(polygons[0]), expected, Math.abs(expected) * 1e-6);
});

it("clipConvex: clips a subject polygon crossing the antimeridian", () => {
  const viewport = rad([[170, -30], [190, -30], [190, -60], [170, -60]]);
  const clip = geoClipConvex(viewport);
  const subject = polygon(rad([[160, -40], [-160, -40], [-160, -50], [160, -50]]));
  const {polygons} = clipCollect(subject, clip);
  assert.strictEqual(polygons.length, 1);
  const expected = geoArea(polygon([[170, -40], [-170, -40], [-170, -50], [170, -50]]));
  assertInDelta(areaOf(polygons[0]), expected, Math.abs(expected) * 1e-6);
});

it("clipConvex: clips a line crossing the antimeridian without breaking it", () => {
  const viewport = rad([[170, -30], [190, -30], [190, -60], [170, -60]]);
  const clip = geoClipConvex(viewport);
  const line = {type: "LineString", coordinates: rad([[160, -45], [-160, -45]])};
  const {lines} = clipCollect(line, clip);
  assert.strictEqual(lines.length, 1);
  assert.strictEqual(lines[0].length, 2);
  assertInDelta(lines[0][0], [170 * radians, -45 * radians], 1e-6);
  assertInDelta(lines[0][1], [190 * radians, -45 * radians], 1e-6);
});

it("clipConvex: handles a clip polygon containing the pole", () => {
  // A viewport reaching the south pole, with continuous longitudes.
  const viewport = rad([[170, -90], [190, -90], [190, -30], [170, -30]]);
  const clip = geoClipConvex(viewport);
  const subject = polygon(rad([[160, -80], [-160, -80], [-160, -90], [160, -90]]));
  const {polygons} = clipCollect(subject, clip);
  assert.strictEqual(polygons.length, 1);
  const expected = geoArea(polygon([[170, -80], [-170, -80], [-170, -90], [170, -90]]));
  assertInDelta(areaOf(polygons[0]), expected, Math.abs(expected) * 1e-6);
});

it("clipConvex: clips a line near the pole without breaking it", () => {
  const viewport = rad([[170, -90], [190, -90], [190, -30], [170, -30]]);
  const clip = geoClipConvex(viewport);
  const line = {type: "LineString", coordinates: rad([[160, -85], [-160, -85]])};
  const {lines} = clipCollect(line, clip);
  assert.strictEqual(lines.length, 1);
  assert.strictEqual(lines[0].length, 2);
  assertInDelta(lines[0][0], [170 * radians, -85 * radians], 1e-6);
  assertInDelta(lines[0][1], [190 * radians, -85 * radians], 1e-6);
});

it("clipConvex: handles MultiPolygon input", () => {
  const clip = geoClipConvex(rad(triangle));
  const multi = {
    type: "MultiPolygon",
    coordinates: [
      [rad([[-10, -10], [10, -10], [0, 10], [-10, -10]])],
      [rad([[-10, -80], [10, -80], [0, -70], [-10, -80]])]
    ]
  };
  const {polygons} = clipCollect(multi, clip);
  assert.strictEqual(polygons.length, 1);
  assert.deepStrictEqual(polygons[0][0], rad([[-10, -10], [10, -10], [0, 10]]));
});

it("clipConvex: drops sliver rings with area < 1e-12", () => {
  const clip = geoClipConvex(rad([[-1, -1], [1, -1], [0, 1]]));
  const tiny = polygon(rad([[0, 0], [1e-6, 0], [0, 1e-6]])); // area ~5e-13
  const {polygons} = clipCollect(tiny, clip);
  assert.strictEqual(polygons.length, 0);
  const small = polygon(rad([[0, 0], [1e-4, 0], [0, 1e-4]])); // area ~5e-9
  const {polygons: kept} = clipCollect(small, clip);
  assert.strictEqual(kept.length, 1);
});

it("clipConvex: throws on degenerate input", () => {
  const message = /convex polygon requires at least 3 distinct vertices/;
  assert.throws(() => geoClipConvex([]), message);
  assert.throws(() => geoClipConvex(rad([[0, 0]])), message);
  assert.throws(() => geoClipConvex(rad([[0, 0], [1, 1]])), message);
  // Duplicate (coincident) vertices.
  assert.throws(() => geoClipConvex(rad([[0, 0], [0, 0], [1, 1]])), message);
  assert.throws(() => geoClipConvex(rad([[0, 0], [1, 0], [1, 1], [0, 0]])), message);
  // Collinear vertices.
  assert.throws(() => geoClipConvex(rad([[0, 0], [1, 0], [2, 0]])), message);
  // Self-intersecting (bowtie).
  assert.throws(() => geoClipConvex(rad([[0, 0], [1, 1], [1, 0], [0, 1]])), message);
  // Non-convex.
  assert.throws(() => geoClipConvex(rad([[0, 0], [2, 0], [1, 0.5], [2, 1], [0, 1]])), message);
  // Longitude span of a full turn.
  assert.throws(() => geoClipConvex(rad([[-180, -80], [180, -80], [180, 80], [-180, 80]])), message);
});

it("clipConvex: geoPath area matches the spherical area", () => {
  const small = [[0, 0], [0.2, 0], [0, 0.2]];
  const projection = geoAzimuthalEqualArea()
      .scale(1)
      .translate([0, 0])
      .precision(1e-3)
      .preclip(geoClipConvex(rad(small)));
  const path = geoPath(projection);
  const expected = Math.abs(sphericalTriangleArea(rad(small)));
  assertInDelta(Math.abs(path.area(worldRectangle)), expected, 1e-6);
});

