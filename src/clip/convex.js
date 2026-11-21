import {merge} from "d3-array";
import {abs, atan2, cos, epsilon, pi, quarterPi, sin, tau} from "../math.js";
import clipBuffer from "./buffer.js";
import clipRejoin from "./rejoin.js";

var errorMessage = "convex polygon requires at least 3 distinct vertices";

// Clips a stream to a convex polygon specified as a list of [lambda, phi]
// points (spherical coordinates in radians, not closed). Clipping is done in
// the equirectangular (lambda, phi) plane, generalizing geoClipRectangle:
// edges are straight lines in lambda and phi. The polygon may cross the
// antimeridian if given with continuous longitudes (e.g. 170° to 190°), and
// may touch a pole. Both winding orders are accepted; the output preserves
// the input winding.
export default function(points) {
  if (!points || points.length < 3) throw new Error(errorMessage);

  var n = points.length,
      v = [],
      i, j, x, y, z;

  for (i = 0; i < n; ++i) v.push([+points[i][0], +points[i][1]]);

  // Edges, and the polygon’s longitude extent (must be less than a full turn).
  var ex = new Array(n), ey = new Array(n), el2 = new Array(n), el = new Array(n),
      nx = new Array(n), ny = new Array(n),
      lambda0 = Infinity, lambda1 = -Infinity, w = 0;

  for (i = 0; i < n; ++i) {
    j = (i + 1) % n;
    x = v[j][0] - v[i][0];
    y = v[j][1] - v[i][1];
    if ((el2[i] = x * x + y * y) < 1e-24) throw new Error(errorMessage); // coincident vertices
    el[i] = Math.sqrt(el2[i]);
    ex[i] = x, ey[i] = y;
    if (v[i][0] < lambda0) lambda0 = v[i][0];
    if (v[i][0] > lambda1) lambda1 = v[i][0];
  }
  if (lambda1 - lambda0 >= tau - epsilon) throw new Error(errorMessage);

  // Strict convexity: every consecutive pair of edges must turn the same way.
  for (i = 0; i < n; ++i) {
    j = (i + 1) % n;
    z = ex[i] * ey[j] - ey[i] * ex[j];
    if (abs(z) < 1e-12) throw new Error(errorMessage); // collinear or duplicate
    if (!w) w = z > 0 ? 1 : -1;
    else if (z * w < 0) throw new Error(errorMessage); // non-convex or self-intersecting
  }

  var lc = (lambda0 + lambda1) / 2; // mid-longitude: center of the unwrapped frame

  // Outward edge normals: a point p is inside iff n_i . (p - v_i) <= 0.
  for (i = 0; i < n; ++i) nx[i] = w * ey[i], ny[i] = -w * ex[i];

  // Shifts a longitude into the polygon’s frame, [lc - pi, lc + pi).
  function shiftLambda(x) {
    return x + tau * Math.ceil((lc - x) / tau - 0.5);
  }

  function visible(x, y) {
    x = shiftLambda(x);
    for (var i = 0; i < n; ++i) {
      if (nx[i] * (x - v[i][0]) + ny[i] * (y - v[i][1]) > 1e-12) return false;
    }
    return true;
  }

  // Cyrus–Beck clipping of the segment a -> b against the polygon; a and b
  // are mutated to the clipped endpoints. Returns true if any part is visible.
  function clipSegment(a, b) {
    var dx = b[0] - a[0], dy = b[1] - a[1],
        ax = a[0], ay = a[1],
        t0 = 0, t1 = 1, i, s, q, t;
    for (i = 0; i < n; ++i) {
      s = nx[i] * dx + ny[i] * dy;
      q = nx[i] * (ax - v[i][0]) + ny[i] * (ay - v[i][1]);
      if (s > 0) { // exiting
        t = -q / s;
        if (t < t0) return false;
        if (t < t1) t1 = t;
      } else if (s < 0) { // entering
        t = -q / s;
        if (t > t1) return false;
        if (t > t0) t0 = t;
      } else if (q > 0) { // parallel and outside
        return false;
      }
    }
    if (t0 > t1) return false;
    if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
    if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
    return true;
  }

  // Parameterizes the polygon boundary in input order: vertex i has
  // parameter i, and points on edge i interpolate linearly towards i + 1.
  function boundaryParam(p) {
    var x = shiftLambda(p[0]), y = p[1],
        best = 0, bestd = Infinity, i, d, s;
    for (i = 0; i < n; ++i) {
      d = abs(nx[i] * (x - v[i][0]) + ny[i] * (y - v[i][1])) / el[i];
      if (d < bestd) bestd = d, best = i;
    }
    s = ((x - v[best][0]) * ex[best] + (y - v[best][1]) * ey[best]) / el2[best];
    return best + (s < 0 ? 0 : s > 1 ? 1 : s);
  }

  // Emits the polygon vertices between two boundary points, walking the
  // boundary in the input order (direction = 1) or in reverse (-1).
  function interpolate(from, to, direction, stream) {
    var i, a, b, k;
    if (from == null) {
      if (direction > 0) { for (i = 0; i < n; ++i) stream.point(v[i][0], v[i][1]); }
      else { for (i = n - 1; i >= 0; --i) stream.point(v[i][0], v[i][1]); }
      return;
    }
    a = boundaryParam(from);
    b = boundaryParam(to);
    if (direction > 0) {
      if (b <= a) b += n;
      for (k = Math.floor(a) + 1; k < Math.ceil(b); ++k) stream.point(v[k % n][0], v[k % n][1]);
    } else {
      if (b >= a) b -= n;
      for (k = Math.ceil(a) - 1; k > Math.floor(b); --k) {
        i = k % n;
        if (i < 0) i += n;
        stream.point(v[i][0], v[i][1]);
      }
    }
  }

  function compareIntersection(a, b) {
    return boundaryParam(a.x) - boundaryParam(b.x);
  }

  function clip(stream) {
    var activeStream = stream,
        bufferStream = clipBuffer(),
        segments,
        polygon,
        ring,
        x__, y__, // first point, as received
        x_, y_, v_, // previous point, in the polygon’s frame
        xRaw_, // previous point’s longitude, as received
        first,
        clean,
        allVisible; // whether every ring of the polygon was fully visible

    var clipStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: polygonStart,
      polygonEnd: polygonEnd,
      sphere: sphere
    };

    function point(x, y) {
      x = shiftLambda(x);
      if (visible(x, y)) activeStream.point(x, y);
    }

    // Winding number of the subject polygon around the clip polygon’s first
    // vertex (the reference point where the boundary parameterization wraps,
    // matching the convention of clipRejoin). Each ring’s longitude span
    // realigns the reference’s 2π branch first, so that rings covering the
    // whole longitude range (e.g. a world rectangle) are handled correctly.
    function polygonInside() {
      var winding = 0, i, j, m, r, a, b, lo, hi, kx, ky = v[0][1];
      for (i = 0; i < polygon.length; ++i) {
        if ((m = (r = polygon[i]).length) < 2) continue;
        lo = hi = r[0][0];
        for (j = 1; j < m; ++j) {
          if (r[j][0] < lo) lo = r[j][0];
          if (r[j][0] > hi) hi = r[j][0];
        }
        kx = v[0][0] + tau * Math.round(((lo + hi) / 2 - v[0][0]) / tau);
        // If the reference lands exactly on the ring’s seam, nudge it inwards.
        if (abs(kx - lo) < 1e-12) kx = lo + (hi - lo) * 1e-9;
        else if (abs(kx - hi) < 1e-12) kx = hi - (hi - lo) * 1e-9;
        for (j = 1; j < m; ++j) {
          a = r[j - 1], b = r[j];
          if (a[1] <= ky) {
            if (b[1] > ky && (b[0] - a[0]) * (ky - a[1]) > (kx - a[0]) * (b[1] - a[1])) ++winding;
          } else {
            if (b[1] <= ky && (b[0] - a[0]) * (ky - a[1]) < (kx - a[0]) * (b[1] - a[1])) --winding;
          }
        }
      }
      return winding;
    }

    function polygonStart() {
      activeStream = bufferStream;
      segments = [];
      polygon = [];
      clean = true;
      allVisible = true;
    }

    function polygonEnd() {
      var startInside = polygonInside(),
          visible = (segments = merge(segments)).length,
          cleanInside = clean && startInside && !allVisible;
      if (cleanInside || visible) {
        stream.polygonStart();
        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }
        if (visible) {
          clipRejoin(segments, compareIntersection, startInside, interpolate, stream);
        }
        stream.polygonEnd();
      }
      activeStream = stream;
      segments = polygon = ring = null;
      allVisible = false;
    }

    function lineStart() {
      clipStream.point = linePoint;
      if (polygon) polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = xRaw_ = NaN;
    }

    function lineEnd() {
      if (segments) {
        linePoint(x__, y__); // close the ring
        clipRing(ring);
        segments.push(bufferStream.result());
      }
      clipStream.point = point;
      if (v_) activeStream.lineEnd();
    }

    // Computes the longitude of the next point in the polygon’s frame,
    // keeping edges continuous: edges spanning the full 2π (e.g. the top and
    // bottom of a world rectangle) are preserved as full sweeps, other edges
    // take the short way around.
    function frameLambda(x) {
      var rdx;
      if (first) return shiftLambda(x);
      rdx = x - xRaw_;
      if (abs(rdx) > pi && abs(abs(rdx) - tau) > epsilon) return x + tau * Math.round((x_ - x) / tau);
      return x_ + rdx;
    }

    function linePoint(x, y) {
      var xs = frameLambda(x);
      if (polygon) {
        ring.push([xs, y]);
      } else {
        linePointVisible(xs, y);
      }
      if (first) {
        x__ = x, y__ = y;
        first = false;
      }
      x_ = xs, y_ = y, xRaw_ = x;
    }

    // Incremental clipping for lines (non-polygon mode).
    function linePointVisible(xs, y) {
      var v = visible(xs, y);
      if (first) {
        if (v) {
          activeStream.lineStart();
          activeStream.point(xs, y);
        }
      } else {
        clipEdge([x_, y_], [xs, y], v_, v);
      }
      v_ = v;
    }

    // Emits the visible part of the edge a -> b, where va and vb are the
    // visibilities of its endpoints.
    function clipEdge(a, b, va, vb) {
      if (va && vb) {
        activeStream.point(b[0], b[1]);
      } else {
        var A = [a[0], a[1]], B = [b[0], b[1]],
            k = tau * Math.round((lc - (A[0] + B[0]) / 2) / tau);
        A[0] += k, B[0] += k;
        if (clipSegment(A, B)) {
          A[0] -= k, B[0] -= k;
          if (!va) {
            activeStream.lineStart();
            activeStream.point(A[0], A[1]);
          }
          activeStream.point(B[0], B[1]);
          if (!vb) activeStream.lineEnd();
          clean = false;
        } else if (vb) {
          activeStream.lineStart();
          activeStream.point(b[0], b[1]);
          clean = false;
        }
      }
    }

    // Clips a buffered ring (polygon mode). Rings that span the whole
    // longitude range wrap the cylinder; their seam edges (pairs of
    // coincident but reversed edges, such as the left and right sides of a
    // world rectangle) are not real boundaries and are dropped.
    function clipRing(r) {
      var m = r.length, i, j, a, b, va, vb, skip,
          lo = Infinity, hi = -Infinity,
          haveA = false, firstVa, lastVa;
      if (m < 2) return;
      for (i = 0; i < m; ++i) {
        if (r[i][0] < lo) lo = r[i][0];
        if (r[i][0] > hi) hi = r[i][0];
      }
      if (hi - lo >= tau - epsilon) {
        skip = new Array(m);
        for (i = 1; i < m; ++i) {
          for (j = i + 1; j < m; ++j) {
            if (congruent(r[i - 1], r[j]) && congruent(r[i], r[j - 1])) {
              skip[i] = skip[j] = true;
              break;
            }
          }
        }
      }
      for (i = 1; i < m; ++i) {
        if (skip && skip[i]) {
          haveA = false;
          continue;
        }
        b = r[i];
        vb = visible(b[0], b[1]);
        if (!vb) allVisible = false;
        if (!haveA) {
          a = r[i - 1];
          va = visible(a[0], a[1]);
          if (!va) allVisible = false;
          if (va) {
            activeStream.lineStart();
            activeStream.point(a[0], a[1]);
          }
          if (firstVa === undefined) firstVa = va;
          haveA = true;
        }
        clipEdge(a, b, va, vb);
        a = b, va = vb;
        lastVa = vb;
      }
      if (firstVa && lastVa) bufferStream.rejoin();
    }

    function sphere() {
      stream.polygonStart();
      stream.lineStart();
      interpolate(null, null, 1, stream);
      stream.lineEnd();
      stream.polygonEnd();
    }

    return clipStream;
  }

  // Drops sliver rings (area < 1e-12) produced by clipping.
  function filter(stream) {
    var rings, ring;
    return {
      point: function(x, y) { if (ring) ring.push([x, y]); else stream.point(x, y); },
      lineStart: function() { if (rings) ring = []; else stream.lineStart(); },
      lineEnd: function() {
        if (rings) rings.push(ring), ring = null;
        else stream.lineEnd();
      },
      polygonStart: function() { rings = []; },
      polygonEnd: function() {
        var kept = null, i, j, m;
        for (i = 0; i < rings.length; ++i) {
          if (abs(ringArea(rings[i])) >= 1e-12) (kept || (kept = [])).push(rings[i]);
        }
        if (kept) {
          stream.polygonStart();
          for (i = 0; i < kept.length; ++i) {
            stream.lineStart();
            for (j = 0, m = kept[i].length; j < m; ++j) stream.point(kept[i][j][0], kept[i][j][1]);
            stream.lineEnd();
          }
          stream.polygonEnd();
        }
        rings = null;
      },
      sphere: function() { stream.sphere(); }
    };
  }

  return function(stream) {
    return clip(filter(stream));
  };
}

// Two points are congruent if they coincide modulo a full turn in longitude.
function congruent(a, b) {
  var d = abs(a[0] - b[0]) % tau;
  return (d < epsilon || d > tau - epsilon) && abs(a[1] - b[1]) < epsilon;
}

// Signed spherical area of a ring of [lambda, phi] points in radians,
// using the same spherical-excess formula as geoArea.
function ringArea(ring) {
  var m = ring.length, i = 0, p,
      lambda0, cosPhi0, sinPhi0,
      lambda, phi, cosPhi, sinPhi,
      dLambda, sdLambda, adLambda, k, u,
      sum = 0;
  if (!m) return 0;
  p = ring[0];
  lambda0 = p[0], phi = p[1] / 2 + quarterPi, cosPhi0 = cos(phi), sinPhi0 = sin(phi);
  while (++i <= m) {
    p = ring[i % m];
    lambda = p[0], phi = p[1] / 2 + quarterPi, cosPhi = cos(phi), sinPhi = sin(phi);
    dLambda = lambda - lambda0;
    sdLambda = dLambda >= 0 ? 1 : -1;
    adLambda = sdLambda * dLambda;
    k = sinPhi0 * sinPhi;
    u = cosPhi0 * cosPhi + k * cos(adLambda);
    sum += atan2(k * sdLambda * sin(adLambda), u);
    lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
  }
  return 2 * sum;
}
