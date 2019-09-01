# robust-predicates

Fast robust predicates for computational geometry in JavaScript (`orient2d`, `orient3d`, `incircle`, `insphere`). Provides reliable 2D and 3D point orientation tests that are not susceptible to floating point errors without sacrificing performance. A modern port of [Jonathan R Shewchuk's C code](https://www.cs.cmu.edu/~quake/robust.html), an industry standard since 1996.

<a href="https://observablehq.com/@mourner/non-robust-arithmetic-as-art"><img height="200" src="predicates.png" /></a>

_Figure: non-robust vs robust `orient2d` test for points within a tiny range (2<sup>-42</sup>)._

[![Build Status](https://travis-ci.com/mourner/robust-predicates.svg?branch=master)](https://travis-ci.com/mourner/robust-predicates)
[![Simply Awesome](https://img.shields.io/badge/simply-awesome-brightgreen.svg)](https://github.com/mourner/projects)
[![Browser Build](https://badgen.net/bundlephobia/minzip/robust-predicates)](https://unpkg.com/robust-predicates)

## [Demo](https://observablehq.com/@mourner/non-robust-arithmetic-as-art)

## API

### `orient2d(ax, ay, bx, by, cx, cy)`

Return a positive value if the points `a`, `b`, and `c` occur in counterclockwise order;
a negative value if they occur in clockwise order; and zero if they are collinear.
The result is also a rough approximation of twice the signed area of the triangle defined by the three points.

### `incircle(ax, ay, bx, by, cx, cy, dx, dy)`

Return a positive value if the point `d` lies inside the circle passing through `a`, `b`, and `c`;
a negative value if it lies outside; and zero if the four points are cocircular.
The points `a`, `b`, and `c` must be in counterclockwise order, or the sign of the result will be reversed.

### `orient3d(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)`

Return a positive value if the point `d` lies below the plane passing through `a`, `b`, and `c`;
"below" is defined so that `a`, `b`, and `c` appear in counterclockwise order when viewed from above the plane.
Returns a negative value if `d` lies above the plane. Returns zero if the points are coplanar.
The result is also a rough approximation of six times the signed volume of the tetrahedron defined by the four points.

### `insphere(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez)`

Return a positive value if the point `e` lies inside the sphere passing through `a`, `b`, `c`, and `d`;
a negative value if it lies outside; and zero if the five points are cospherical.
The points `a`, `b`, `c`, and `d` must be ordered so that they have a positive orientation
(as defined by `orient3d`), or the sign of the result will be reversed.

### `orient2dfast`, `orient3dfast`, `incirclefast`, `inspherefast`

Simple, approximate, non-robust versions of predicates above. Use when robustness isn't needed.

## Example

```js
import {orient2d} from 'robust-predicates';

const ccw = orient2d(ax, ay, bx, by, cx, cy) > 0;
````

## Install

Install with `npm install robust-predicates` or `yarn add robust-predicates`, or use one of the browser builds:

- all predicates: [predicates.min.js](https://unpkg.com/robust-predicates@1.1.0/umd/predicates.min.js)
- `orient2d`, `orient2dfast`: [orient2d.min.js](https://unpkg.com/robust-predicates@1.1.0/umd/orient2d.min.js)
- `orient3d`, `orient3dfast`: [orient3d.min.js](https://unpkg.com/robust-predicates@1.1.0/umd/orient3d.min.js)
- `incircle`, `incirclefast`: [incircle.min.js](https://unpkg.com/robust-predicates@1.1.0/umd/incircle.min.js)
- `insphere`, `inspherefast`: [insphere.min.js](https://unpkg.com/robust-predicates@1.1.0/umd/insphere.min.js)

## License

Since the original code by J. Shewchuk is in the public domain, this port follows the same choice. See [Unlicense](https://unlicense.org)
