# robust-predicates

Fast robust predicates for computational geometry in JavaScript. A modern port of [Jonathan R Shewchuk's C code](https://www.cs.cmu.edu/~quake/robust.html), which has been an industry standard since 1996. Uses code generation to produce highly efficient JavaScript code. _A work in progress:_

- [x] `orient2d`
- [ ] `orient3d`
- [x] `incircle`
- [ ] `insphere`

[![Build Status](https://travis-ci.com/mourner/robust-predicates.svg?branch=master)](https://travis-ci.com/mourner/robust-predicates)
[![Simply Awesome](https://img.shields.io/badge/simply-awesome-brightgreen.svg)](https://github.com/mourner/projects)

## API

### orient2d(ax, ay, bx, by, cx, cy)

Return a positive value if the points `a`, `b`, and `c` occur in counterclockwise order;
a negative value if they occur in clockwise order; and zero if they are collinear.
The result is also a rough approximation of twice the signed area of the triangle defined by the three points.

### incircle(ax, ay, bx, by, cx, cy, dx, dy)

Return a positive value if the point `d` lies inside the circle passing through `a`, `b`, and `c`;
a negative value if it lies outside; and zero if the four points are cocircular.
The points `a`, `b`, and `c` must be in counterclockwise order, or the sign of the result will be reversed.

## Example

```js
import {orient2d} from 'robust-predicates';

const ccw = orient2d(ax, ay, bx, by, cx, cy) > 0;
````

## Install

Install with `npm install robust-predicates` or `yarn add robust-predicates`, or use one of the browser builds:

- [predicates.min.js](https://unpkg.com/robust-predicates@0.1.0/umd/predicates.min.js) (all predicates)
- [orient2d.min.js](https://unpkg.com/robust-predicates@0.1.0/umd/orient2d.min.js) (`orient2d`)
- [incircle.min.js](https://unpkg.com/robust-predicates@0.1.0/umd/incircle.min.js) (`incircle`)

## License

Since the original code by J. Shewchuk is in the public domain, this port follows the same choice. See [Unlicense](https://unlicense.org)
